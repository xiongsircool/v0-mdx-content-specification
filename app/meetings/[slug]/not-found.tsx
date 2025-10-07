import Link from "next/link"
import { Card, CardContent } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Calendar } from "lucide-react"

export default function NotFound() {
  return (
    <div className="min-h-screen bg-background flex items-center justify-center">
      <Card className="max-w-md w-full mx-4">
        <CardContent className="p-8 text-center">
          <Calendar className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
          <h1 className="text-2xl font-bold mb-2">会议未找到</h1>
          <p className="text-muted-foreground mb-6">
            抱歉，您要查找的会议不存在或已被移除。
          </p>
          <div className="space-y-2">
            <Button asChild className="w-full">
              <Link href="/meetings">返回会议列表</Link>
            </Button>
            <Button variant="outline" asChild className="w-full">
              <Link href="/">返回首页</Link>
            </Button>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}